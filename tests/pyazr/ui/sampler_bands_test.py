#!/usr/bin/env python3
"""The sampler and the bands as jobs of the web UI (docs/upgrade/web_ui.md section 11.7), on a copy of T3 with two
free parameters (the level 2- 19.0 MeV: energy_5 and width_5_1; every other level and width fixed; there the
posterior is Gaussian: its standard deviations on a grid equal the Laplace ones to 1e-4), at the environment's
threads (one under ctest):

  1. emcee from the best point of the engine's MIGRAD: the chain's standard deviation of each parameter over the
     fit's (the Laplace approximation, HESSE) within 10 %; the chain in output/chain_<id>.npy (stored steps,
     walkers, parameters) with its JSON summary; the page receives quantiles, split R-hat and the autocorrelation
     time, never the draws;
  2. a stop between steps ends the sampler with the chain so far; a sampler started from it continues its chain
     (the old steps first, unchanged);
  3. the linear band from the fit's covariance and the band from 400 draws of the chain agree: the median ratio of
     their half widths over every point within 10 %; every data set and extrapolation has a band on the points of
     the Results payload; a band job without a covariance or without a source job is refused (422);
  4. the start walkers lie inside the fit's limits; a start outside them is refused naming the parameter.

The numbers are printed (measured on the product).
Run: python3 tests/pyazr/ui/sampler_bands_test.py
"""
import json
import os
import sys
import time

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uitest  # noqa: E402
from uitest import copy_project, make_client, wait_job, wait_until  # noqa: E402

FREE_LEVEL = "2- 19.0"  # T3's level whose energy and first width stay free


def two_parameter_t3(dest):
    """A copy of T3 with the energy and the first width of FREE_LEVEL free, everything else fixed."""
    import pyazr
    path = copy_project("T3", dest)
    p = pyazr.Project.open(path)
    for level in p.levels():
        if level.label != FREE_LEVEL:
            p.set_level(level.label, fixed=True)
        for k, c in enumerate(level.channels):
            if level.label != FREE_LEVEL or k > 0:
                p.set_channel(level.label, c.pair, c.l, c.s, fixed=True)
    p.save(path)
    return path


@pytest.fixture(scope="module")
def t3(tmp_path_factory):
    path = two_parameter_t3(tmp_path_factory.mktemp("sampler"))
    with make_client(path) as client:
        client.path = path
        yield client


def start(client, **body):
    r = client.post("/api/jobs", json=body)
    assert r.status_code == 200, r.text
    return r.json()


@pytest.fixture(scope="module")
def fit_and_chain(t3):
    fit = wait_job(t3, start(t3, kind="fit", method="engine_migrad", x0="text")["id"])
    assert fit["state"] == "converged" and len(fit["result"]["errors"]) == 2, fit["result"]
    t0 = time.perf_counter()
    sample = wait_job(t3, start(t3, kind="sample", method="emcee", x0={"job": fit["id"]},
                                settings={"steps": 2500, "walkers": 8, "spread": 1e-4})["id"])
    print(f"\n  sampler: {time.perf_counter() - t0:.1f} s, {sample['result']['message']}")
    return fit, sample


def test_chain_sd_against_laplace_and_what_the_page_receives(t3, fit_and_chain):
    fit, sample = fit_and_chain
    r = sample["result"]
    assert sample["state"] == "finished", r
    for name, laplace in fit["result"]["errors"].items():
        ratio = r["errors"][name] / laplace
        print(f"  {name}: chain sd {r['errors'][name]:.6g}, Laplace sd {laplace:.6g}, ratio {ratio:.4f}, "
              f"tau {r['tau'][name]:.3g}, split R-hat {r['rhat'][name]:.4f}, effective draws "
              f"{r['effective_draws'][name]:.0f}")
        assert abs(ratio - 1) < 0.10, (name, ratio)
        assert r["quantiles"][name]["16"] < r["quantiles"][name]["50"] < r["quantiles"][name]["84"]
        assert abs(r["rhat"][name] - 1) < 0.1
    chain = np.load(r["chain"])
    assert chain.shape == (2500, 8, 2) and r["chain"].endswith(f"chain_{sample['id']}.npy")
    assert os.path.dirname(r["chain"]) == os.path.join(os.path.dirname(t3.path), "output")
    with open(r["summary"]) as f:
        summary = json.load(f)
    assert summary["names"] == list(fit["result"]["errors"]) and summary["discard"] == 1250
    # the parametrization of the sampled session, as a prior file's #parametrization tag (pyazr.priors.BASES)
    import pyazr
    from pyazr.priors import BASES, chain_basis
    with pyazr.Session(t3.path) as s:
        assert summary["parametrization"] == s.basis and summary["parametrization"] in BASES
    assert chain_basis(r["chain"]) == summary["parametrization"]
    flat = chain[1250:].reshape(-1, 2)
    assert np.allclose([r["mean"][n] for n in summary["names"]], flat.mean(axis=0), rtol=0, atol=0)
    assert all(not isinstance(v, list) or len(v) < 10 for v in r.values()), "the page receives no draws"
    assert sample["progress"]["step"] == 2500 and sample["progress"]["tau"] > 0
    assert 0 < sample["progress"]["acceptance"] < 1


def test_stop_between_steps_and_resume(t3, fit_and_chain):
    fit, _ = fit_and_chain
    job = start(t3, kind="sample", method="emcee", x0={"job": fit["id"]}, settings={"steps": 100000})
    path = f"/api/jobs/{job['id']}"
    assert wait_until(lambda: ((t3.get(path).json().get("progress") or {}).get("step") or 0) >= 5, 120)
    t0 = time.perf_counter()
    t3.post(f"{path}/stop")
    stopped = wait_job(t3, job["id"])
    print(f"\n  sampler stop: request to 'stopped' {time.perf_counter() - t0:.3f} s, "
          f"{stopped['result']['reports_after_stop']} step(s) after it")
    assert stopped["state"] == "stopped" and stopped["result"]["reports_after_stop"] <= 1
    old = np.load(stopped["result"]["chain"])
    assert old.shape[0] >= 5
    more = wait_job(t3, start(t3, kind="sample", method="emcee", x0={"job": job["id"]},
                              settings={"steps": 20})["id"])
    new = np.load(more["result"]["chain"])
    assert more["result"]["continues"] == job["id"]
    assert new.shape == (old.shape[0] + 20,) + old.shape[1:]
    assert np.array_equal(new[:old.shape[0]], old)


def test_killed_sampler_continues_from_its_saved_chain(t3, fit_and_chain):
    """A sampler killed (or interrupted by a restart: no result either) continues from its chain file as last saved
    (every SAVE_EVERY stored steps); the page offers the sampler beside a sampler job (its follows)."""
    from pyazr.ui.jobs import SAVE_EVERY, job_row
    fit, _ = fit_and_chain
    job = start(t3, kind="sample", method="emcee", x0={"job": fit["id"]}, settings={"steps": 100000})
    path = f"/api/jobs/{job['id']}"
    assert wait_until(lambda: ((t3.get(path).json().get("progress") or {}).get("step") or 0) > SAVE_EVERY + 2, 300)
    t3.post(f"{path}/kill")
    killed = wait_job(t3, job["id"])
    assert killed["state"] == "killed" and not killed.get("result"), killed
    saved = os.path.join(os.path.dirname(t3.path), "output", f"chain_{job['id']}.npy")
    old = np.load(saved)
    assert old.shape[0] >= SAVE_EVERY and old.shape[0] % SAVE_EVERY == 0
    more = wait_job(t3, start(t3, kind="sample", method="emcee", x0={"job": job["id"]}, settings={"steps": 10})["id"])
    new = np.load(more["result"]["chain"])
    assert more["result"]["continues"] == job["id"] and more["result"]["continued_from"] == saved
    assert new.shape == (old.shape[0] + 10,) + old.shape[1:] and np.array_equal(new[:old.shape[0]], old)
    assert ("sample", "emcee") in job_row("sample", "emcee").follows


def test_linear_band_and_band_from_draws_agree(t3, fit_and_chain):
    fit, sample = fit_and_chain
    linear = wait_job(t3, start(t3, kind="bands", method="covariance", x0={"job": fit["id"]})["id"])
    draws = wait_job(t3, start(t3, kind="bands", method="draws", x0={"job": sample["id"]},
                               settings={"draws": 400})["id"])
    assert linear["state"] == draws["state"] == "finished", (linear["result"], draws["result"])
    # which parameters each band varied, from its result: the level parameters (linear), every free one (draws)
    free = list(fit["result"]["errors"])
    assert draws["result"]["varied"] == draws["result"]["free"] == free
    assert linear["result"]["free"] == free and set(linear["result"]["varied"]) <= set(free)
    ev = t3.get("/api/evaluation").json()
    names = [d["name"] for d in ev["datasets"]] + [d["name"] for d in ev["extrapolations"]]
    points = {d["name"]: len(d["e_cm"]) for d in ev["datasets"] + ev["extrapolations"]}
    ratios = []
    for a, b in zip(linear["result"]["bands"]["sets"], draws["result"]["bands"]["sets"]):
        assert a["name"] == b["name"] and len(a["lower"]) == len(b["upper"]) == points[a["name"]]
        half_linear = (np.array(a["upper"]) - np.array(a["lower"])) / 2
        half_draws = (np.array(b["upper"]) - np.array(b["lower"])) / 2
        ratios.append(half_draws[half_linear > 0] / half_linear[half_linear > 0])
    assert [s["name"] for s in linear["result"]["bands"]["sets"]] == names
    ratios = np.concatenate(ratios)
    median = float(np.median(ratios))
    print(f"\n  bands: draws over linear half width at {ratios.size} points: median {median:.4f}, 5-95 % "
          f"{np.percentile(ratios, 5):.3f}-{np.percentile(ratios, 95):.3f}, within 10 % at "
          f"{np.mean(np.abs(ratios - 1) < 0.1):.1%} of the points")
    assert abs(median - 1) < 0.10
    listed = next(j for j in t3.get("/api/jobs").json() if j["id"] == linear["id"])
    assert listed["result"]["bands"] == "GET /api/jobs/{id}", "the list leaves the arrays out"
    scipy = wait_job(t3, start(t3, kind="fit", method="scipy_lbfgsb", x0="text")["id"])
    r = t3.post("/api/jobs", json={"kind": "bands", "method": "covariance", "x0": {"job": scipy["id"]}})
    assert r.status_code == 422 and "covariance" in r.json()["message"], r.text
    r = t3.post("/api/jobs", json={"kind": "bands", "method": "draws"})
    assert r.status_code == 422 and "chain" in r.json()["message"], r.text


def test_start_walkers_inside_the_limits_and_a_start_outside_refused():
    from pyazr.errors import ParameterError
    from pyazr.ui.jobs import _start_walkers
    rng = np.random.default_rng(1)
    bounds = [(0.0, 1.0), (-np.inf, np.inf)]
    p0 = _start_walkers(lambda x: 0.0, np.array([0.999, 5.0]), bounds, ["a", "b"], 50, 0.01, rng)
    assert p0.shape == (50, 2) and np.all(p0[:, 0] <= 1.0) and np.all(p0[:, 0] >= 0.0)
    with pytest.raises(ParameterError, match=r"a = 1\.5 \(limits \[0, 1\]\)"):
        _start_walkers(lambda x: 0.0, np.array([1.5, 5.0]), bounds, ["a", "b"], 4, 0.01, rng)


if __name__ == "__main__":
    uitest.run(__file__)
