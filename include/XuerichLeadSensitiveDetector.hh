#ifndef __XENON10PLEADSENSITIVEDETECTOR_H__
#define __XENON10PLEADSENSITIVEDETECTOR_H__

#include <G4VSensitiveDetector.hh>

#include "XuerichLeadHit.hh"

class G4Step;
class G4HCofThisEvent;

class XuerichLeadSensitiveDetector: public G4VSensitiveDetector
{
public:
	XuerichLeadSensitiveDetector(G4String hName);
	~XuerichLeadSensitiveDetector();

	void Initialize(G4HCofThisEvent *pHitsCollectionOfThisEvent);
	G4bool ProcessHits(G4Step *pStep, G4TouchableHistory *pHistory);
	void EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent);

private:
	XuerichLeadHitsCollection* m_pLeadHitsCollection;

	map<int,G4String> m_hParticleTypes;
};

#endif
