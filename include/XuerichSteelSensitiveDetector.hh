#ifndef __XENON10PSTEELSENSITIVEDETECTOR_H__
#define __XENON10PSTEELSENSITIVEDETECTOR_H__

#include <G4VSensitiveDetector.hh>

#include "XuerichSteelHit.hh"

class G4Step;
class G4HCofThisEvent;

class XuerichSteelSensitiveDetector: public G4VSensitiveDetector
{
public:
	XuerichSteelSensitiveDetector(G4String hName);
	~XuerichSteelSensitiveDetector();

	void Initialize(G4HCofThisEvent *pHitsCollectionOfThisEvent);
	G4bool ProcessHits(G4Step *pStep, G4TouchableHistory *pHistory);
	void EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent);

private:
	XuerichSteelHitsCollection* m_pSteelHitsCollection;

	map<int,G4String> m_hParticleTypes;
};

#endif
