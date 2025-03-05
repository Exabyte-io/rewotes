import React from 'react';
import { render } from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';

it('should render a heading', () => {
  expect(true).toBe(true);
});

it('should render a heading', () => {
  const { getByText } = render(<h1>Hello, world!</h1>);
  expect(getByText(/Hello, world!/)).toBeInTheDocument();
});
