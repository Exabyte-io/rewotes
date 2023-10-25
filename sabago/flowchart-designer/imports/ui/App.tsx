import React, { useState } from 'react';
import FlowchartDesigner from './Designer';
import InitialConditionsForm from './InitialCondForm';

const App = () => {
  const [dopingLevel, setDopingLevel] = useState<number>(0); // Let's say this is in atoms/cm^3 for simplicity.
  const [disableNDoping, setDisableNDoping] = useState<boolean>(false);
  const [disablePDoping, setDisablePDoping] = useState<boolean>(false);
  const [desiredCond, setDesiredCond] = useState<number>(0);
  const [formSubmitted, setFormSubmitted] = useState<boolean>(false);

  const handleFormSubmit = () => {
    setFormSubmitted(true);
  };
  return (
    <div style={{ width: '100%', height: '100vh' }}>
      <h2 style={{ color: '#00165A' }}>Welcome to MAT3RA!</h2>
      <h3 style={{ color: '#00228D' }}>
        With this feature, you can create a flowchart for doping Si wafers to the desired conductivity
      </h3>
      <p style={{ color: '#00228D' }}>
        See `flowchartdesigner.md` (at the root of this project on Github) for more details
      </p>
      <h4 style={{ color: '#00228D' }}>
        To start, set the initial doping amount and the desired conductivity in the form below:
      </h4>
      <InitialConditionsForm
        formSubmitted={formSubmitted}
        initialDopingLevel={dopingLevel}
        desiredConductivity={desiredCond}
        onDopingChange={setDopingLevel}
        disableNDoping={setDisableNDoping}
        disablePDoping={setDisablePDoping}
        onConductivityChange={setDesiredCond}
        onSubmit={() => {
          handleFormSubmit();
        }}
      />
      <FlowchartDesigner
        desiredConductivity={desiredCond}
        setDesiredConductivity={setDesiredCond}
        isFormSubmitted={formSubmitted}
        setIsFormSubmitted={setFormSubmitted}
        dopingLevel={dopingLevel}
        setDopingLevel={setDopingLevel}
        disabledNDoping={disableNDoping}
        disabledPDoping={disablePDoping}
        setDisableNdoping={setDisableNDoping}
        setDisablePDoping={setDisablePDoping}
      />
    </div>
  );
};

export default App;
