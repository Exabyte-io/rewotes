import React from 'react';
import { SaveButton } from './styles';

export interface ICFormProps {
  initialDopingLevel: number;
  onDopingChange: React.Dispatch<React.SetStateAction<number>>;
  desiredConductivity: number;
  onConductivityChange: React.Dispatch<React.SetStateAction<number>>;
  disablePDoping: (disable: boolean) => void;
  disableNDoping: (disable: boolean) => void;
  onSubmit: () => void;
  formSubmitted: boolean;
}

const InitialConditionsForm = ({
  initialDopingLevel,
  desiredConductivity,
  onDopingChange,
  onConductivityChange,
  onSubmit,
  disablePDoping,
  disableNDoping,
  formSubmitted,
}: ICFormProps) => {
  const handleDopingChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    onDopingChange(Number(event.target.value));
  };

  const handlePDopingChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    onDopingChange(Number(event.target.value));
    disableNDoping(true);
  };

  const handleNDopingChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    onDopingChange(Number(event.target.value));
    disablePDoping(true);
  };

  const handleConductivityChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    // const range = event.target.value.split('-').map(Number);
    onConductivityChange(Number(event.target.value));
  };
  return (
    <form
      style={{ marginBottom: '20px', display: 'flex' }}
      onSubmit={(e) => {
        e.preventDefault();
        onSubmit();
      }}
    >
      <fieldset>
        <legend>Initial Doping Amount:</legend>
        <label>
          <input
            type="radio"
            value="0"
            checked={initialDopingLevel === 0}
            name="doping"
            onChange={handleDopingChange}
            disabled={formSubmitted}
          />
          0 atoms/cm^3 (Undoped Si)
        </label>
        <label>
          <input
            type="radio"
            value="10"
            checked={initialDopingLevel === 10}
            name="doping"
            onChange={handlePDopingChange}
            disabled={formSubmitted}
          />
          10 atoms/cm^3 (p-type doped)
        </label>
        <label>
          <input
            type="radio"
            value="30"
            checked={initialDopingLevel === 30}
            name="doping"
            onChange={handleNDopingChange}
            disabled={formSubmitted}
          />
          30 atoms/cm^3 (n-type doped)
        </label>
      </fieldset>

      {initialDopingLevel === 0 ? (
        <fieldset>
          <legend>Desired Conductivity:</legend>
          <label>
            <input
              type="radio"
              value="0.4"
              checked={desiredConductivity === 0.4}
              name="conductivity"
              onChange={handleConductivityChange}
              disabled={formSubmitted}
            />
            0.4 S/m
          </label>
          <label>
            <input
              type="radio"
              value="0.82"
              checked={desiredConductivity === 0.82}
              name="conductivity"
              onChange={handleConductivityChange}
              disabled={formSubmitted}
            />
            0.82 S/m
          </label>
        </fieldset>
      ) : initialDopingLevel === 10 ? (
        <fieldset>
          <legend>Desired Conductivity:</legend>
          <label>
            <input
              type="radio"
              value="0.05"
              checked={desiredConductivity === 0.05}
              name="conductivity"
              onChange={handleConductivityChange}
              disabled={formSubmitted}
            />
            0.05 S/m
          </label>
          <label>
            <input
              type="radio"
              value="0.4"
              checked={desiredConductivity === 0.4}
              name="conductivity"
              onChange={handleConductivityChange}
              disabled={formSubmitted}
            />
            0.4 S/m
          </label>
        </fieldset>
      ) : (
        <fieldset>
          <legend>Desired Conductivity:</legend>
          <label>
            <input
              type="radio"
              value="0.2"
              checked={desiredConductivity === 0.2}
              name="conductivity"
              onChange={handleConductivityChange}
              disabled={formSubmitted}
            />
            0.2 S/m
          </label>
          <label>
            <input
              type="radio"
              value="1.7"
              checked={desiredConductivity === 1.7}
              name="conductivity"
              onChange={handleConductivityChange}
            />
            1.7 S/m
          </label>
        </fieldset>
      )}
      <SaveButton type="submit" disabled={desiredConductivity === 0 || formSubmitted}>
        {' '}
        Save values
      </SaveButton>
    </form>
  );
};

export default InitialConditionsForm;
