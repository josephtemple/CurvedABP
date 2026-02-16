#pragma once
#include "state.h"
#include "vec2.h"

constexpr double PI = 3.14159265358979323846;

Vec2 gradV(double x, double y, double A);
void step(Simulation& sim);
