#include "cellMechanics.h"

CellMechanics::CellMechanics(HemoCellField & cellfield) : cellConstants(CommonCellConstants::CommonCellConstantsConstructor(cellfield)) {}
