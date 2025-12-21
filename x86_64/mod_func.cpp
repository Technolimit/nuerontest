#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern "C" void _kdr_reg(void);
extern "C" void _naf_reg(void);
extern "C" void _vecst_reg(void);

extern "C" void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"kdr.mod\"");
    fprintf(stderr, " \"naf.mod\"");
    fprintf(stderr, " \"vecst.mod\"");
    fprintf(stderr, "\n");
  }
  _kdr_reg();
  _naf_reg();
  _vecst_reg();
}
