#ifndef __XENON10PSCINTSENSITIVEDETECTOR_H__
#define __XENON10PSCINTSENSITIVEDETECTOR_H__

#include <G4VSensitiveDetector.hh>

#include "XuerichNaIHit.hh"

class G4Step;
class G4HCofThisEvent;

class XuerichNaISensitiveDetector: public G4VSensitiveDetector
{
public:
	XuerichNaISensitiveDetector(G4String hName);
	~XuerichNaISensitiveDetector();

	void Initialize(G4HCofThisEvent *pHitsCollectionOfThisEvent);
	G4bool ProcessHits(G4Step *pStep, G4TouchableHistory *pHistory);
	void EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent);

private:
	XuerichNaIHitsCollection* m_pNaIHitsCollection;

	map<int,G4String> m_hParticleTypes;
};

#endif // __XENON10PLXESENSITIVEDETECTOR_H__
