import React from 'react';
import { act, render, screen } from '@testing-library/react';
// import '@testing-library/jest-dom/extend-expect';
import App from '../../imports/ui/App';
import InitialConditionsForm from '../../imports/ui/InitialCondForm';
import FlowchartDesigner from '../../imports/ui/Designer';

jest.mock('../../imports/ui/InitialCondForm', () => {
  const mock = jest.fn();
  (mock as jest.Mock).mockReturnValue(null); // You can use mockImplementation if you need more complex behavior
  return mock;
});

jest.mock('../../imports/ui/Designer', () => {
  const mock = jest.fn();
  (mock as jest.Mock).mockReturnValue(null);
  return mock;
});

describe('App Component', () => {
  beforeEach(() => {
    jest.clearAllMocks();
  });

  it('renders without crashing', () => {
    render(<App />);
    expect(screen.getByText('Welcome to MAT3RA!')).toBeInTheDocument();
    expect(
      screen.getByText(
        'With this feature, you can create a flowchart for doping Si wafers to the desired conductivity',
      ),
    ).toBeInTheDocument();
  });

  it('renders InitialConditionsForm', () => {
    render(<App />);
    expect(InitialConditionsForm).toHaveBeenCalled();
  });

  it('renders FlowchartDesigner', () => {
    render(<App />);
    expect(FlowchartDesigner).toHaveBeenCalled();
  });

  it('updates formSubmitted state on form submission', () => {
    render(<App />);

    // Check the isFormSubmitted prop is passed to FlowchartDesigner
    const getFormSubmittedProp = () => {
      const designerCalls = (FlowchartDesigner as jest.Mock).mock.calls;
      const lastCallProps = designerCalls.length > 0 ? designerCalls[designerCalls.length - 1][0] : undefined;
      return lastCallProps ? lastCallProps.isFormSubmitted : undefined;
    };

    // Initially, form should not be submitted
    expect(getFormSubmittedProp()).toBe(false);

    act(() => {
      const submitFunction = (InitialConditionsForm as jest.Mock).mock.calls[0][0].onSubmit;
      submitFunction();
    });

    // Check if form has been submitted
    expect(getFormSubmittedProp()).toBe(true);
  });
});
