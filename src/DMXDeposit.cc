#include "DMXDeposit.hh"

Deposits::DMXDeposit::DMXDeposit()
{
}

Deposits::DMXDeposit::~DMXDeposit()
{
}

void Deposits::DMXDeposit::PlaceDeposit(G4Event* event, std::unique_ptr<const Deposits::EnergyDeposit> deposit)
{
    fDeposits.emplace_back(std::make_pair(event->GetEventID(), std::move(deposit)));
}
