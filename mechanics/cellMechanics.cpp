#include "cellMechanics.h"

CellMechanics::CellMechanics(HemoCellField & cellfield, Config & modelCfg_) : cellConstants(CommonCellConstants::CommonCellConstantsConstructor(cellfield, modelCfg_)) {}
