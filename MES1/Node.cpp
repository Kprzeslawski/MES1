#include "Node.h"

void Node::set_bc(){
	BC = 1;
}


double Node::get_x() { return this->x; }
double Node::get_y() { return this->y; }
int Node::bc() { return this->BC; }