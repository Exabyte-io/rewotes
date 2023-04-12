import { createContext, useContext } from 'react';

export interface Settings {
  theme?: string;
  //TODO: add atomColor
  // other user settings 
}

interface SettingsContextType {
  settings: Settings;
  updateSettings: (settings: Settings) => void;
}

export const SettingsContext = createContext<SettingsContextType | undefined>(undefined);

export function useSettings(): SettingsContextType {
    const context = useContext(SettingsContext);
  
    if (!context) {
      throw new Error('useSettings must be used within a SettingsProvider');
    }
  
    return context;
  }