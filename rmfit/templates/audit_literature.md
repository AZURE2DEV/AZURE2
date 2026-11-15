# Literature review for the $reaction R-matrix fit (`$azr`)

$date.  Requested by rmfit for the campaign `$campaign_dir`.  Rules and report format at the end.

Compound system and channels in the fit (from the model's pairs):

$pairs

## What to find

1. **Earlier evaluations of the same compound system**: what data they used, which channels, radii, energy limits,
   how they treated normalizations, backgrounds and energy shifts; where their parameters differ from this fit's
   level scheme (table of J^pi, E_x, widths with citations).
2. **Measurements the fit does not yet use**: channels open in the window without data in the segment table; newer
   or independent measurements of the sets that are in; absolute measurements that could anchor a free scale.
3. **Method papers** on the specific difficulty the campaign names (readme below): normalization treatment, Park vs
   Brune, background poles, data weighting, resolution/target convolution.
4. **Level scheme**: ENSDF (`pyazr.nds.fetch_levels`) vs the model's levels in the window: missing levels with a
   known J^pi, J^pi conflicts, energies off by more than the quoted uncertainty.

Tag every source by access: [full] (text read), [abstract], [cited] (through another paper; name it), [snippet].
Mark your own reasoning "Inference:".  Nothing from memory.  Paywalled papers: `~/.claude/skills/rclone-box-drive`
(Box `Analysis/<reaction>/`), NASA ADS for pre-1995 scans.  Worked example: an archived campaign.

The output is a SHORT list of actionable ideas, each with a cheap test in this campaign (a dataset-off run, a
frozen evaluation with a changed radius, a level add at a literature energy), not a survey.  Use the same report
table: `Where` = the citation, `Evidence` = what the paper measured or did, `Cheapest check` = the test here.

## The model's segment table (what the fit uses)

$segment_table

## Already known (do not re-report)

$known

## The evaluator's framing (campaign readme, first lines)

$readme_head

$common
