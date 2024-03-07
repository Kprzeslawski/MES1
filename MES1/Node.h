#pragma once
class Node
{
	float x, y;
	float t;
	int BC;

public:

	Node(float x, float y, float t, int bc) :x(x), y(y), t(t), BC(bc) {}
	Node(float x, float y, float t) :x(x), y(y), t(t), BC(0) {}
	Node() {};

	void set_bc();
	double get_x();
	double get_y();
	int bc();

};

