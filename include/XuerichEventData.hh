#ifndef __XENON10PEVENTDATA_H__
#define __XENON10PEVENTDATA_H__

#include <string>
#include <vector>

using std::string;
using std::vector;

class XuerichEventData
{
public:
	 XuerichEventData();
	~XuerichEventData();

public:
	void Clear();

public:
	int m_iEventId;								// the event ID

	int m_iNbSteps;								// number of energy depositing steps
	int m_iNbStepsNaI;						
	int m_iNbStepsTeflon;						
	int m_iNbStepsDeadXe;						
	int m_iNbStepsGXe;							
	int m_iNbStepsSteel;							
	int m_iNbStepsAl;							
	int m_iNbStepsLead;							
	int m_iNbStepsRest;							

	int m_iNbTopPmtHits;						// number of top pmt hits
	int m_iNbBottomPmtHits;						// number of bottom pmt hits
	vector<int> *m_pPmtHits;					// number of photon hits per pmt

	vector<int> *m_pTrackId;					// id of the particle
	vector<int> *m_pParentId;					// id of the parent particle

	vector<string> *m_pParticleType;			// type of particle
	vector<string> *m_pParticleTypeNaI;			// type of particle

	vector<string> *m_pParentType;				// type of particle
	vector<string> *m_pCreatorProcess;			// interaction
	vector<string> *m_pDepositingProcess;		// energy depositing process

	vector<float> *m_pX;						// position of the step in Xe target
	vector<float> *m_pY;
	vector<float> *m_pZ;

	vector<float> *m_pXNaI;						// position of the step in NaIillator
	vector<float> *m_pYNaI;
	vector<float> *m_pZNaI;

	vector<float> *m_pKineticEnergy;			// particle kinetic energy after the step			
    vector<double>* m_pPreKinNrg;
	vector<double>* m_pPostKinNrg;

	vector<float> *m_pEnergyDeposited; 			// energy deposited in the step
	vector<float> *m_pEnergyDepositedNaI;
	vector<float> *m_pEnergyDepositedTeflon;
	vector<float> *m_pEnergyDepositedDeadXe;
	vector<float> *m_pEnergyDepositedGXe;
	vector<float> *m_pEnergyDepositedSteel;
	vector<float> *m_pEnergyDepositedAl;
	vector<float> *m_pEnergyDepositedLead;
	vector<float> *m_pEnergyDepositedRest;

	//float m_fTotalEnergyDepositedER;			// total energy deposited in the NaISD
	//float m_fTotalEnergyDepositedNR;
	//float m_fTotalEnergyDepositedNaIER;
	//float m_fTotalEnergyDepositedNaINR;
	float m_fTotalEnergyDeposited;
	float m_fTotalEnergyDepositedNaI;
	float m_fTotalEnergyDepositedTeflon;
	float m_fTotalEnergyDepositedDeadXe;
	float m_fTotalEnergyDepositedGXe;
	float m_fTotalEnergyDepositedSteel;
	float m_fTotalEnergyDepositedAl;
	float m_fTotalEnergyDepositedLead;
	float m_fTotalEnergyDepositedRest;

	vector<float> *m_pTime;						// time of the step
	vector<float> *m_pTimeNaI;				
	vector<float> *m_pTimeTeflon;				
	vector<float> *m_pTimeDeadXe;
	vector<float> *m_pTimeGXe;
	vector<float> *m_pTimeSteel;
	vector<float> *m_pTimeAl;
	vector<float> *m_pTimeLead;
	vector<float> *m_pTimeRest;	

	vector<string> *m_pPrimaryParticleType;		// type of particle
	float m_fPrimaryX;							// position of the primary particle
	float m_fPrimaryY;
	float m_fPrimaryZ;	
};

#endif // __XENON10PEVENTDATA_H__

