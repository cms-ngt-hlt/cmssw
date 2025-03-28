// -*- C++ -*-
#ifndef FWCore_Framework_DependentRecordImplementation_h
#define FWCore_Framework_DependentRecordImplementation_h
//
// Package:     Framework
// Class  :     DependentRecordImplementation
//
/**\class edm::eventsetup::DependentRecordImplementation

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Author:      Chris Jones
// Created:     Fri Apr 29 10:03:54 EDT 2005
//

// system include files
#include <sstream>
#include <type_traits>

// user include files
#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/EventSetupImpl.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/NoRecordException.h"
#include "FWCore/Framework/interface/DependentRecordTag.h"
#include "FWCore/Utilities/interface/mplVector.h"

//This is here only because too many modules depend no
// getting this header from this file (before EventSetupImpl)
#include "FWCore/Framework/interface/EventSetup.h"

// forward declarations
namespace edm {
  namespace eventsetup {

    template <class RecordT, class ListT>
    class DependentRecordImplementation : public EventSetupRecordImplementation<RecordT>, public DependentRecordTag {
    public:
      DependentRecordImplementation() {}
      using list_type = ListT;
      //virtual ~DependentRecordImplementation();

      // ---------- const member functions ---------------------
      template <class DepRecordT>
      const DepRecordT getRecord() const {
        //Make sure that DepRecordT is a type in ListT
        static_assert(
            (list_type::template contains<DepRecordT>()),
            "Trying to get a Record from another Record where the second Record is not dependent on the first Record.");
        try {
          EventSetup const eventSetupT{
              this->eventSetup(), this->transitionID(), this->getTokenIndices(), *this->esParentContext()};
          return eventSetupT.get<DepRecordT>();
        } catch (cms::Exception& e) {
          std::ostringstream sstrm;
          sstrm << "While getting dependent Record from Record " << this->key().type().name();
          e.addContext(sstrm.str());
          throw;
        }
      }

      template <class DepRecordT>
      std::optional<DepRecordT> tryToGetRecord() const {
        //Make sure that DepRecordT is a type in ListT
        static_assert(
            (list_type::template contains<DepRecordT>()),
            "Trying to get a Record from another Record where the second Record is not dependent on the first Record.");
        EventSetup const eventSetupT{
            this->eventSetup(), this->transitionID(), this->getTokenIndices(), *this->esParentContext()};
        return eventSetupT.tryToGet<DepRecordT>();
      }

      using EventSetupRecordImplementation<RecordT>::getHandle;

      template <typename ProductT, typename DepRecordT>
      ESHandle<ProductT> getHandle(ESGetToken<ProductT, DepRecordT> const& iToken) const {
        //Make sure that DepRecordT is a type in ListT
        static_assert((list_type::template contains<DepRecordT>()),
                      "Trying to get a product with an ESGetToken specifying a Record from another Record where the "
                      "second Record is not dependent on the first Record.");
        return getRecord<DepRecordT>().getHandle(iToken);
      }

      using EventSetupRecordImplementation<RecordT>::getTransientHandle;

      template <typename ProductT, typename DepRecordT>
      ESTransientHandle<ProductT> getTransientHandle(ESGetToken<ProductT, DepRecordT> const& iToken) const {
        //Make sure that DepRecordT is a type in ListT
        static_assert((list_type::template contains<DepRecordT>()),
                      "Trying to get a product with an ESGetToken specifying a Record from another Record where the "
                      "second Record is not dependent on the first Record.");
        return getRecord<DepRecordT>().getTransientHandle(iToken);
      }

      using EventSetupRecordImplementation<RecordT>::get;

      template <typename ProductT, typename DepRecordT>
      ProductT const& get(ESGetToken<ProductT, DepRecordT> const& iToken) const {
        //Make sure that DepRecordT is a type in ListT
        static_assert((list_type::template contains<DepRecordT>()),
                      "Trying to get a product with an ESGetToken specifying a Record from another Record where the "
                      "second Record is not dependent on the first Record.");
        return getRecord<DepRecordT>().get(iToken);
      }

      template <typename ProductT, typename DepRecordT>
      ProductT const& get(ESGetToken<ProductT, DepRecordT>& iToken) const {
        return get(const_cast<ESGetToken<ProductT, DepRecordT> const&>(iToken));
      }

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------

    private:
      // ---------- member data --------------------------------
    };

  }  // namespace eventsetup
}  // namespace edm

#endif
