#ifndef LOADBALANCER_CPP
#define LOADBALANCER_CPP

#include "loadBalancer.h"

LoadBalancer::LoadBalancer(HemoCell & hemocell_) : hemocell(hemocell_) {

}

double LoadBalancer::calculateFractionalLoadImbalance() {
  return 0.;
}

void LoadBalancer::doLoadBalance() {

}
#endif
