//____________________________________________________________________________
/*!

\class   genie::MCModel

\brief

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#include <sstream>

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "ReWeight/MCModel.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;

//___________________________________________________________________________
MCModel::MCModel()
{
  fName = "unamed mc model";
  this->Initialize();
}
//___________________________________________________________________________
MCModel::MCModel(string name)
{
  fName = name;
  this->Initialize();
}
//___________________________________________________________________________
MCModel::MCModel(const MCModel & model)
{
  this->Copy(model);
}
//___________________________________________________________________________
MCModel::~MCModel()
{

}
//___________________________________________________________________________
void MCModel::Copy(const MCModel & model)
{

}
//___________________________________________________________________________
void MCModel::UseXSecAlg(const ProcessInfo & proc, const AlgId & algid)
{
  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * alg = 
              dynamic_cast<const XSecAlgorithmI *> 
                         (algf->GetAlgorithm(algid.Name(), algid.Config()));

  string key = this->BuildKey(proc);

  fXSecModelList.insert(
                 map<string, const XSecAlgorithmI *>::value_type(key, alg));
}
//___________________________________________________________________________
void MCModel::UseXSecAlg(
    const ProcessInfo & proc, const InitialState & init, const AlgId & algid)
{
  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * alg = 
              dynamic_cast<const XSecAlgorithmI *> 
                         (algf->GetAlgorithm(algid.Name(), algid.Config()));

  string key = this->BuildKey(proc, init);

  fXSecModelList.insert(
                 map<string, const XSecAlgorithmI *>::value_type(key, alg));
}
//___________________________________________________________________________
const XSecAlgorithmI * MCModel::XSecAlg(const Interaction * interaction) const
{
  const ProcessInfo &  proc = interaction->GetProcessInfo();
  const InitialState & init = interaction->GetInitialState();

  string key = this->BuildKey(proc, init);

  if(fXSecModelList.count(key) == 1) 
  {
    map<string, const XSecAlgorithmI *>::const_iterator iter;
    iter = fXSecModelList.find(key);

    const XSecAlgorithmI * alg = iter->second;
    return alg;
  }

  key = this->BuildKey(proc);

  if(fXSecModelList.count(key) == 1) 
  {
    map<string, const XSecAlgorithmI *>::const_iterator iter;
    iter = fXSecModelList.find(key);

    const XSecAlgorithmI * alg = iter->second;
    return alg;
  }

  LOG("ReWeight", pWARN) 
           << "No cross section model for the input interaction";
  return 0;
}
//___________________________________________________________________________
void MCModel::Initialize(void)
{
  fXSecModelList.clear();
}
//___________________________________________________________________________
string MCModel::BuildKey(const ProcessInfo & proc) const
{
  ostringstream key;
  key << "PROC:" << proc.AsString();
  return key.str();
}
//___________________________________________________________________________
string MCModel::BuildKey(
                   const ProcessInfo & proc, const InitialState & init) const
{
  ostringstream key;
  key << "PROC:" << proc.AsString() << ";INIT:"  <<  init.AsString();
  return key.str();
}
//___________________________________________________________________________

