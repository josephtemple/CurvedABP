// dynamics.h
//
// gradient function and simulation stepping.

#pragma once
#include "state.h"
#include "vec.h"

Vec2 gradV(double x, double y, double A);
void step(Simulation& sim);
