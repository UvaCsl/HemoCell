#ifndef LOADBALANCER_H
#define LOADBALANCER_H

class LoadBalancer;

#include "hemocell_internal.h"
#include "hemocell.h"

class LoadBalancer {
  public:
  LoadBalancer(HemoCell & hemocell);
  double calculateFractionalLoadImbalance();
  void doLoadBalance();

  private:
  HemoCell & hemocell;
};

#endif
