// dynamics.h
//
// gradient function and simulation stepping.

#pragma once
#include "state.h"
#include "vec.h"

constexpr double PI = 3.14159265358979323846;

Vec2 gradV(double x, double y, double A);
void step(Simulation& sim);
