#ifndef WASAUNPACKHELPER
#define WASAUNPACKHELPER

#include "WASAUnpackBranch.hh"


class EventWASAUnpack
{
public:
 S4TQ s4tq;
 S4MWDC s4mwdc;
 S4WFD s4wfd;
 S2TQ1 s2tq1;
 S2MDC s2mdc;
 S2WFD123 s2wfd123;
 S2TQ2 s2tq2;
 S2WFD45 s2wfd45;
 S2Fiber s2fiber;
 S2CsI s2csi;
 FRSTPC frstpc;

 EventWASAUnpack();
 ~EventWASAUnpack() = default;
}