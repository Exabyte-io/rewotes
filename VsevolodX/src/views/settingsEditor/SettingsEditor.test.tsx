import React from 'react';
import { render, fireEvent, screen } from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import { SettingsProvider } from '../../context/SettingsContext';
import SettingsEditor from './SettingsEditor';

const renderSettingsEditor = () => {
  return render(
    <SettingsProvider>
      <SettingsEditor />
    </SettingsProvider>
  );
};

describe('SettingsEditor', () => {
  it('toggles the theme when the button is clicked', () => {
    renderSettingsEditor();
    const themeButtonStateDark = screen.getByRole('button', { name: /Dark Theme/i });
    expect(themeButtonStateDark).toBeInTheDocument();
    expect(screen.getByRole('heading', { level: 4 })).toHaveTextContent('Settings');
    fireEvent.click(themeButtonStateDark);
    expect(themeButtonStateDark).toHaveTextContent('Light Theme');

    const themeButtonStateLight = screen.getByRole('button', { name: /Light Theme/i });
    expect(themeButtonStateLight).toBeInTheDocument();
    expect(screen.getByRole('heading', { level: 4 })).toHaveTextContent('Settings');
    fireEvent.click(themeButtonStateLight);
    expect(themeButtonStateLight).toHaveTextContent('Dark Theme');

  });

  it('loads atoms preferences when the button is clicked', async () => {
    global.fetch = jest.fn(() =>
    Promise.resolve({
        ok: true,
        json: () => Promise.resolve({ atomsDisplayData: [{ element: 'H', color: '#ffffff' }] }),
    } as Response)
);

    renderSettingsEditor();
    const loadButton = screen.getByRole('button', { name: /Load atoms preferences/i });
    expect(loadButton).toBeInTheDocument();
    fireEvent.click(loadButton);
    expect(global.fetch).toHaveBeenCalledTimes(1);
    await screen.findByText('H');
  });

  it('toggles the edit mode when the button is clicked', () => {
  renderSettingsEditor();
  //Button by default lets enable 3D mode
  const editButtonState1 = screen.getByRole('button', { name: /Edit 3D/i });
  expect(editButtonState1).toBeInTheDocument();

  fireEvent.click(editButtonState1);
  expect(editButtonState1).toHaveTextContent('Edit Source');

  //Button changes state to toggle to Source Mode
  const editButtonState2 = screen.getByRole('button', { name: /Edit Source/i });
  expect(editButtonState2).toBeInTheDocument();

  fireEvent.click(editButtonState2);
  expect(editButtonState2).toHaveTextContent('Edit 3D');
});
});
