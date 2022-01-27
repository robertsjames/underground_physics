#ifndef DMXDeposit_h
#define DMXDeposit_h 1

#include "G4RunManager.hh"
#include "NEST.hh"

namespace Deposits
{
struct EnergyDeposit
{
    const G4double fEdep;
    const G4ThreeVector fPos;
    const G4double fTime;
    const NEST::INTERACTION_TYPE fType;

    EnergyDeposit(const G4double edep, const G4ThreeVector pos,
                  const G4double time, const NEST::INTERACTION_TYPE type)
      : fEdep(edep), fPos(pos), fTime(time), fType(type) {}
};

class DMXDeposit
{
public:
    DMXDeposit();
    ~DMXDeposit();

    void PlaceDeposit(G4Event* event, std::unique_ptr<const Deposits::EnergyDeposit> deposit);

private:
    std::vector<std::pair<G4int, std::unique_ptr<const Deposits::EnergyDeposit>>> fDeposits;
};
}

#endif
