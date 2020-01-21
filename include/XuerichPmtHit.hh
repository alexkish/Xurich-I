#ifndef __XENON10PPMTHIT_H__
#define __XENON10PPMTHIT_H__

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4Allocator.hh>
#include <G4ThreeVector.hh>

class XuerichPmtHit: public G4VHit
{
public:
	XuerichPmtHit();
	~XuerichPmtHit();
	XuerichPmtHit(const XuerichPmtHit &);
	const XuerichPmtHit & operator=(const XuerichPmtHit &);
	G4int operator==(const XuerichPmtHit &) const;

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	void Draw();
	void Print();

public:
	void SetPosition(G4ThreeVector hPosition) { m_hPosition = hPosition; }
	void SetTime(G4double dTime) { m_dTime = dTime; }
	void SetPmtNb(G4int iPmtNb) { m_iPmtNb = iPmtNb; }

	G4ThreeVector GetPosition() { return m_hPosition; }
	G4double GetTime() { return m_dTime; }
	G4int GetPmtNb() { return m_iPmtNb; }

private:
	G4ThreeVector m_hPosition;
	G4double m_dTime;
	G4int m_iPmtNb;
};

typedef G4THitsCollection<XuerichPmtHit> XuerichPmtHitsCollection;

extern G4Allocator<XuerichPmtHit> XuerichPmtHitAllocator;

inline void*
XuerichPmtHit::operator new(size_t)
{
	return((void *) XuerichPmtHitAllocator.MallocSingle());
}

inline void
XuerichPmtHit::operator delete(void *pXuerichPmtHit)
{
	XuerichPmtHitAllocator.FreeSingle((XuerichPmtHit*) pXuerichPmtHit);
}

#endif // __XENON10PPMTHIT_H__

