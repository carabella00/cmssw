 // -*- C++ -*-
//
// Package:     SiStripObjects
// Class  :     SiStripDetCabling
// Implementation:
//     <Notes on implementation>
// Original Author:  dkcira
//         Created:  Wed Mar 22 12:24:33 CET 2006
// $Id: SiStripDetCabling.cc,v 1.5 2007/01/29 15:24:35 dkcira Exp $

#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;


SiStripDetCabling::SiStripDetCabling()
{
}


// SiStripDetCabling::SiStripDetCabling(const SiStripDetCabling& rhs)
// {
//    // do actual copying here;
// }


SiStripDetCabling::~SiStripDetCabling()
{
}


SiStripDetCabling::SiStripDetCabling(const SiStripFedCabling& fedcabling) : fullcabling_(), connected_(), detected_(), undetected_()
{
  // --- CONNECTED = have fedid and i2cAddr
  // create fullcabling_, loop over vector of FedChannelConnection, either make new element of map, or add to appropriate vector of existing map element
  // get feds list (vector) from fedcabling object - these are the active FEDs
  const vector<uint16_t>& feds = fedcabling.feds();
  vector<uint16_t>::const_iterator ifed;
  for ( ifed = feds.begin(); ifed != feds.end(); ifed++ ) { // iterate over active feds, get all their FedChannelConnection-s
    const vector<FedChannelConnection>& conns = fedcabling.connections( *ifed );
    vector<FedChannelConnection>::const_iterator iconn;
    for ( iconn = conns.begin(); iconn != conns.end(); iconn++ ) { // loop over FedChannelConnection objects
      addDevices(*iconn, fullcabling_); // leave separate method, in case you will need to add devices also after constructing
      bool have_fed_id = iconn->fedId();
      vector<int> vector_of_connected_apvs;
      if(have_fed_id){ // these apvpairs are seen from the readout
	// there can be at most 6 APVs on one DetId: 0,1,2,3,4,5
        int which_apv_pair = iconn->apvPairNumber(); // APVPair (0,1) for 512 strips and (0,1,2) for 768 strips
	if(iconn->i2cAddr(0)) vector_of_connected_apvs.push_back(2*which_apv_pair + 0); // first apv of the pair
	if(iconn->i2cAddr(1)) vector_of_connected_apvs.push_back(2*which_apv_pair + 1); // second apv of the pair
      }
      if(vector_of_connected_apvs.size() != 0){ // add only is smth. there, obviously
        map<uint32_t, vector<int> > map_of_connected_apvs;
        map_of_connected_apvs.insert(std::make_pair(iconn->detId(),vector_of_connected_apvs));
        addFromSpecificConnection(connected_, map_of_connected_apvs);
      }
    }
  }
  // --- DETECTED = do not have fedid but have i2cAddr
  const std::vector<FedChannelConnection>& detected_fed_connections = fedcabling.detected();
  for(std::vector<FedChannelConnection>::const_iterator idtct = detected_fed_connections.begin(); idtct != detected_fed_connections.end(); idtct++){
    addDevices(*idtct, fullcabling_);
    bool have_fed_id = idtct->fedId();
    vector<int> vector_of_detected_apvs;
    if(! have_fed_id){
      int which_apv_pair = idtct->apvPairNumber(); // APVPair (0,1) for 512 strips and (0,1,2) for 768 strips
      if(idtct->i2cAddr(0)) vector_of_detected_apvs.push_back(2*which_apv_pair + 0); // first apv of the pair
      if(idtct->i2cAddr(1)) vector_of_detected_apvs.push_back(2*which_apv_pair + 1); // second apv of the pair
    }
    if(vector_of_detected_apvs.size() != 0){ // add only is smth. there, obviously
      map<uint32_t,vector<int> > map_of_detected_apvs;
      map_of_detected_apvs.insert(std::make_pair(idtct->detId(),vector_of_detected_apvs));
      addFromSpecificConnection(detected_, map_of_detected_apvs);
    }
  }
  // --- UNDETECTED = have neither fedid nor i2caddr
  const std::vector<FedChannelConnection>& undetected_fed_connections = fedcabling.undetected();
  for(std::vector<FedChannelConnection>::const_iterator iudtct = undetected_fed_connections.begin(); iudtct != undetected_fed_connections.end(); iudtct++){
    addDevices(*iudtct, fullcabling_);
    bool have_fed_id = iudtct->fedId();
    vector<int> vector_of_undetected_apvs;
    if(! have_fed_id){
      int which_apv_pair = iudtct->apvPairNumber(); // APVPair (0,1) for 512 strips and (0,1,2) for 768 strips
      if(iudtct->i2cAddr(0)) vector_of_undetected_apvs.push_back(2*which_apv_pair + 0); // first apv of the pair
      if(iudtct->i2cAddr(1)) vector_of_undetected_apvs.push_back(2*which_apv_pair + 1); // second apv of the pair
    }
    if(vector_of_undetected_apvs.size() != 0){ // add only is smth. there, obviously
      map<uint32_t, vector<int> > map_of_undetected_apvs;
      map_of_undetected_apvs.insert(std::make_pair(iudtct->detId(),vector_of_undetected_apvs));
      addFromSpecificConnection(undetected_, map_of_undetected_apvs);
    }
  }
}


void SiStripDetCabling::addDevices(const FedChannelConnection & conn, std::map< uint32_t, std::vector<FedChannelConnection> >& chosen_connections_ ){ // add to certain connections
// separate method in case need to add devices/cabling after object already constructed
  const uint32_t& fedchannelconnection_detid = conn.detId(); // detId of FedChannelConnection that is looped over
  // is this detid already in the map?
  map< uint32_t, vector<FedChannelConnection> >::iterator detcabl_it = chosen_connections_.find(fedchannelconnection_detid);
  if( ! (detcabl_it==chosen_connections_.end()) ){      // is already in map -- DKwarn : add constistency checks here for example agreement with already existing dcuId?
    ( detcabl_it->second ).push_back(conn); // add one FedChannelConnection to vector mapped to this detid
  }else{                                      // is not in map, create vector with 1 element and put into map
    vector<FedChannelConnection> local_fedchannelconnection; local_fedchannelconnection.push_back(conn);
    chosen_connections_[fedchannelconnection_detid] = local_fedchannelconnection;
  }
}


void SiStripDetCabling::addDevices(const FedChannelConnection & conn){ // by default add to fullcabling_ connections - special case of above class
 addDevices(conn, fullcabling_ ); // add to fullcabling_
}


void SiStripDetCabling::addActiveDetectorsRawIds(std::vector<uint32_t> & vector_to_fill_with_detids ) const{ // replace getActiveDetectorRawIds method - avoid use of static
 for(map< uint32_t, vector<int> >::const_iterator conn_it = connected_.begin(); conn_it!=connected_.end(); conn_it++){
   vector_to_fill_with_detids.push_back(conn_it->first);
 }
 // vector_to_fill_with_detids is empty if no elements in fullcabling_
}


const vector<FedChannelConnection>& SiStripDetCabling::getConnections(uint32_t det_id ) const{ // return all connections corresponding to one det_id
  map< uint32_t, vector<FedChannelConnection> >::const_iterator detcabl_it = fullcabling_.find(det_id); // has to be const_iterator because this function cannot change data members
  if( ! (detcabl_it==fullcabling_.end()) ){  // found detid in fullcabling_
    return ( detcabl_it->second );
  }else{ // DKwarn : is there need for output message here telling det_id does not exist?
    static vector<FedChannelConnection> default_empty_fedchannelconnection;
    return default_empty_fedchannelconnection;
  }
}


const FedChannelConnection& SiStripDetCabling::getConnection( uint32_t det_id, unsigned short apv_pair ) const{
  const vector<FedChannelConnection>& fcconns = getConnections(det_id);
  for(vector<FedChannelConnection>::const_iterator iconn = fcconns.begin(); iconn!=fcconns.end();iconn++){
    if ( (iconn->apvPairNumber()) == apv_pair){ // check if apvPairNumber() of present FedChannelConnection is the same as requested one
      return (*iconn); // if yes, return the FedChannelConnection object
    }
  }
  // if did not match none of the above, return some default value - DKwarn : also output message?
  static FedChannelConnection default_empty_fedchannelconnection;
  return default_empty_fedchannelconnection;
}


const unsigned int SiStripDetCabling::getDcuId( uint32_t det_id ) const{
  const vector<FedChannelConnection>& fcconns = getConnections( det_id );
  if(fcconns.size()!=0) {
     return ( fcconns.at(1) ).dcuId(); // get dcuId of first element - when you build check this consistency
  }
  // default if none of the above is fulfilled
  unsigned int default_zero_value = 0;
  return default_zero_value;
}

    // one can find the nr of apvs from fullcabling_ -> vector<FedChannelConnection> -> size * 2
const uint16_t SiStripDetCabling::nApvPairs(uint32_t det_id) const{
 const vector<FedChannelConnection>& fcconns = getConnections( det_id );
 if(fcconns.size()!=0) {
  return( fcconns.at(1).nApvPairs()); // nr of apvpairs for associated module
 }else{
   return 0;
 }
}


void SiStripDetCabling::addConnected ( std::map<uint32_t, std::vector<int> > & map_to_add_to) const{ // map of detector to list of APVs for APVs seen from FECs and FEDs
  addFromSpecificConnection(map_to_add_to, connected_);
}


void SiStripDetCabling::addDetected( std::map<uint32_t, std::vector<int> > & map_to_add_to) const{ // map of detector to list of APVs for APVs seen neither from FECS or FEDs
  addFromSpecificConnection(map_to_add_to, detected_);
}


void SiStripDetCabling::addUnDetected( std::map<uint32_t, std::vector<int> > & map_to_add_to) const{ // map of detector to list of APVs for APVs seen neither from FECS or FEDs
  addFromSpecificConnection(map_to_add_to, undetected_);
}


void SiStripDetCabling::addNotConnectedAPVs( std::map<uint32_t, std::vector<int> > & map_to_add_to) const{ // map of detector to list of APVs that are not connected - combination of addDetected and addUnDetected
  addFromSpecificConnection(map_to_add_to, detected_);
  addFromSpecificConnection(map_to_add_to, undetected_);
}


void SiStripDetCabling::addFromSpecificConnection( std::map<uint32_t, std::vector<int> > & map_to_add_to, const std::map< uint32_t, vector<int> > & specific_connection) const{
  for(map< uint32_t, vector<int> >::const_iterator conn_it = specific_connection.begin(); conn_it!=specific_connection.end(); conn_it++){
    uint32_t new_detid = conn_it->first;
    vector<int> new_apv_vector = conn_it->second;
    std::map<uint32_t, std::vector<int> >::iterator it = map_to_add_to.find(new_detid);
    if( it == map_to_add_to.end() ){ // detid does not exist in map, add new entry
         std::sort(new_apv_vector.begin(),new_apv_vector.end()); // not very efficient sort, time consuming?
         map_to_add_to.insert(std::make_pair(new_detid,new_apv_vector));
    }else{                    // detid exists already, add to its vector - if its not there already . . .
      vector<int> existing_apv_vector = it->second;
      for(std::vector<int>::iterator inew = new_apv_vector.begin(); inew != new_apv_vector.end(); inew++ ){
        bool there_already = false;
        for(std::vector<int>::iterator iold = existing_apv_vector.begin(); iold != existing_apv_vector.end(); iold++){
           if (*iold == *inew){
              there_already = true;
              break; // leave the loop
           }
        }
        if( ! there_already ){
          existing_apv_vector.push_back(*inew);   
          std::sort(existing_apv_vector.begin(),existing_apv_vector.end()); // not very efficient sort, time consuming?
        }else{
          edm::LogWarning("Logical") << "apv "<<*inew<<" already exists in the detector module "<<new_detid;
        }
      }
    }
  }
}


EVENTSETUP_DATA_REG(SiStripDetCabling);
