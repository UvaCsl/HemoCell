#ifndef BOND_PARTICLE_3D_HH
#define BOND_PARTICLE_3D_HH

#include "BondParticle3D.hh"
using namespace plb;


//namespace trombocit {

/* *************** class BondParticle3D ************************************ */


template<typename T, template<typename U> class Descriptor>
void trombocit::BondParticle3D<T,Descriptor>::advance() {

	bondTime += 1;
	if (particles[0] != NULL and particles[1] != NULL) {
		this->getPosition() = (positions[0] + positions[1]) * 0.5 + dx;
	} else {
		this->getPosition() += dx; //		this->setTag(-1);
	}
	processor = this->getMpiProcessor();
	particles[0] = particles[1] = NULL;
//	pcout << uid << " advanced" << std::endl;
}

template<typename T, template<typename U> class Descriptor>
void trombocit::BondParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Particle3D<T,Descriptor>::serialize(serializer);
    serializer.addValue<plint>(processor);
    serializer.addValue<T>(r);
    serializer.addValue<T>(bondTime);
    serializer.addValues<T,3>(eij);
	serializeString(serializer, uid);
    for (int var = 0; var < 2; ++var) {
		serializer.addValues<T,3>(positions[var]);
//		serializer.addValue<plint>(processors[var]);
		serializer.addValue<plint>(cellId[var]);
		serializer.addValue<plint>(vertexId[var]);
    }
}

template<typename T, template<typename U> class Descriptor>
void trombocit::BondParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    Particle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValue<plint>(processor);
    unserializer.readValue<T>(r);
    unserializer.readValue<T>(bondTime);
    unserializer.readValues<T,3>(eij);
    uid =  unserializeString(unserializer);
    for (int var = 0; var < 2; ++var) {
        unserializer.readValues<T,3>(positions[var]);
//        unserializer.readValue<plint>(processors[var]);
        unserializer.readValue<plint>(cellId[var]);
        unserializer.readValue<plint>(vertexId[var]);
        particles[var] = NULL;
    }
}


template<typename T, template<typename U> class Descriptor>
trombocit::BondParticle3D<T,Descriptor>* trombocit::BondParticle3D<T,Descriptor>::clone() const {
    return new BondParticle3D<T,Descriptor>(*this);
}



//} // namespace trombocit



#endif  // BOND_PARTICLE_3D_HH
