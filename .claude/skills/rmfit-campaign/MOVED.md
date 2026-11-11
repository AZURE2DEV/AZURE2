# rmfit-campaign skill: moved

The `rmfit-campaign` skill (the rmfit process guide for global R-matrix fitting
campaigns) moved out of the AZURE2 repository on 2026-11-07 (DeBoer's decision).
It now lives with the code it describes, in the rmfit repository:

    rmfit/.claude/skills/rmfit-campaign/SKILL.md
    (on the CRC cluster: /groups/rdeboer1/user/rdeboer1/R-matrix/rmfit/.claude/skills/rmfit-campaign/SKILL.md)

and is installed for every session as `~/.claude/skills/rmfit-campaign`.
Edit it there (commit and push in rmfit), not here. The last version in this
repository is `git show 4d16764:.claude/skills/rmfit-campaign/SKILL.md`.

This file is deliberately not named SKILL.md, so that only one skill called
`rmfit-campaign` is ever loaded.
