// Constants
export const BASE_CONDUCTIVITY = 0.072; // S/m, Assumed* for undoped silicon
export const BASE_CARRIER_CONCENTRATION = 1e15; // carriers/cm^3, example value for intrinsic silicon

export const CHARGE_OF_HOLE = 1.6e-19; // Coulombs
export const HOLE_MOBILITY = 450; // cm^2/V·s

export const DESIRED_CONDUCTIVITY = 0.4;
export const MARGIN = 0.05;
export const BATCH_SIZE = 10; // Since this does not affect the simplified conductivity calculation, it is unitless to keep things simple.

export const calculatePCarrierConcentration = (dopingLevel: number) => {
  return BASE_CARRIER_CONCENTRATION + dopingLevel * 1e14;
};

// Calculate Conductivity
export const measurePConductivity = (dopingLevel: number) => {
  const n = calculatePCarrierConcentration(dopingLevel);
  // Calculate conductivity
  const conductivity = n * CHARGE_OF_HOLE * HOLE_MOBILITY;
  return conductivity;
};

export const CHARGE_OF_ELECTRON = 1.6e-19; // Coulombs
export const ELECTRON_MOBILITY = 1350; // cm^2/V·s

export const calculateNCarrierConcentration = (dopingLevel: number) => {
  return BASE_CARRIER_CONCENTRATION + dopingLevel * 1e14;
};

export const measureNConductivity = (dopingLevel: number) => {
  const n = calculateNCarrierConcentration(dopingLevel);
  // Calculate conductivity
  const conductivity = n * CHARGE_OF_ELECTRON * ELECTRON_MOBILITY;
  return conductivity;
};
