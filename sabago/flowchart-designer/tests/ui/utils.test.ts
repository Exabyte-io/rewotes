import {
  BASE_CARRIER_CONCENTRATION,
  CHARGE_OF_HOLE,
  HOLE_MOBILITY,
  calculatePCarrierConcentration,
  measurePConductivity,
  CHARGE_OF_ELECTRON,
  ELECTRON_MOBILITY,
  calculateNCarrierConcentration,
  measureNConductivity,
} from '../../imports/ui/utils';

describe('Silicon doping functions', () => {
  // Tests for P-type carrier concentration
  describe('calculatePCarrierConcentration', () => {
    it('should calculate the P-type carrier concentration for a given doping level', () => {
      const dopingLevel = 2;
      const expectedValue = BASE_CARRIER_CONCENTRATION + dopingLevel * 1e14;
      expect(calculatePCarrierConcentration(dopingLevel)).toBe(expectedValue);
    });
  });

  // Tests for P-type conductivity
  describe('measurePConductivity', () => {
    it('should measure the P-type conductivity for a given doping level', () => {
      const dopingLevel = 2;
      const n = calculatePCarrierConcentration(dopingLevel);
      const expectedValue = n * CHARGE_OF_HOLE * HOLE_MOBILITY;
      expect(measurePConductivity(dopingLevel)).toBe(expectedValue);
    });
  });

  // Tests for N-type carrier concentration
  describe('calculateNCarrierConcentration', () => {
    it('should calculate the N-type carrier concentration for a given doping level', () => {
      const dopingLevel = 2;
      const expectedValue = BASE_CARRIER_CONCENTRATION + dopingLevel * 1e14;
      expect(calculateNCarrierConcentration(dopingLevel)).toBe(expectedValue);
    });
  });

  // Tests for N-type conductivity
  describe('measureNConductivity', () => {
    it('should measure the N-type conductivity for a given doping level', () => {
      const dopingLevel = 2;
      const n = calculateNCarrierConcentration(dopingLevel);
      const expectedValue = n * CHARGE_OF_ELECTRON * ELECTRON_MOBILITY;
      expect(measureNConductivity(dopingLevel)).toBe(expectedValue);
    });
  });
});
