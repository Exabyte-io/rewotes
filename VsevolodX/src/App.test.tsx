import React from 'react';
import { render, screen } from '@testing-library/react';
import App from './App';
import { SettingsProvider } from './context/SettingsContext';

test('renders div with correct className', () => {
  render(
    <SettingsProvider>
      <App />
    </SettingsProvider>
    );

  const divElement = screen.getByTestId('app');
  expect(divElement).toHaveClass('App')
  expect(divElement).toMatch(/bp4-(dark|light)/);
});
