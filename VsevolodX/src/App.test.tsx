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
    
  const divElement = screen.getByRole('main');
  expect(divElement).toHaveClass('App-main')
  //expect(divElement).toMatch(/.*bp4-(dark|light)/);
});
