#ifndef __XENON10PALSENSITIVEDETECTOR_H__
#define __XENON10PALSENSITIVEDETECTOR_H__

#include <G4VSensitiveDetector.hh>

#include "XuerichAlHit.hh"

class G4Step;
class G4HCofThisEvent;

class XuerichAlSensitiveDetector: public G4VSensitiveDetector
{
public:
	XuerichAlSensitiveDetector(G4String hName);
	~XuerichAlSensitiveDetector();

	void Initialize(G4HCofThisEvent *pHitsCollectionOfThisEvent);
	G4bool ProcessHits(G4Step *pStep, G4TouchableHistory *pHistory);
	void EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent);

private:
	XuerichAlHitsCollection* m_pAlHitsCollection;

	map<int,G4String> m_hParticleTypes;
};

#endif
