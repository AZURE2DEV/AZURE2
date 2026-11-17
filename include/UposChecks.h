#ifndef AZURE2_UPOSCHECKS_H
#define AZURE2_UPOSCHECKS_H

class CompoundSystem;
class Config;
class Segment;

/// Whether a data set or extrapolation with a secondary gamma ray (UPOS) can be set up; false after an ERROR in the
/// run's log.  A particle exit: the decay of the heavy partner only.  A photon exit: the capture cascade (the
/// secondary gamma ray of the capture final state after an unobserved primary photon), refused with the RMC
/// formalism and for an analysing power.
bool UPOSAllowed(const CompoundSystem *compound, const Segment &segment, int key, bool isTest, const Config &configure);

/// The WARNING of a UPOS data set given at a lab angle after a particle exit (the angle is converted with the
/// kinematics of the light ejectile, not the gamma ray's); nothing for a capture cascade, whose angle is converted as
/// for any capture gamma ray.
void UPOSLabNote(const CompoundSystem *compound, const Segment &segment, int key, bool isTest, const Config &configure);

#endif
