import * as React from 'react';

export const Canvas = ({ children }: { children: React.ReactNode }) => {
  return <div data-testid="mocked-canvas">{children}</div>;
};

export const useThree = () => {
  return {
    camera: {},
    gl: {
      domElement: document.createElement('canvas'),
    },
  };
};

export const useFrame = () => {};
