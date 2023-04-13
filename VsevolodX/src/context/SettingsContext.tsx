import { createContext, useCallback, useContext, useEffect, useMemo, useState } from 'react';

export interface Settings {
  theme?: string;
  defaultAtomColor?: string;
  atomsData?: AtomData[];
  // other user settings 
}

interface SettingsContextType {
  settings: Settings;
  updateSettings: (settings: Settings) => void;
}

interface AtomData {
  element: string,
  color?: string,
  radius?: string,
}

export const SettingsContext = createContext<SettingsContextType | undefined>(undefined);

export function useSettings(): SettingsContextType {
    const context = useContext(SettingsContext);
  
    if (!context) {
      throw new Error('useSettings must be used within a SettingsProvider');
    }
  
    return context;
  }

  interface SettingsProviderProps {
    children: React.ReactNode;
  }
  

export function SettingsProvider({ children } : SettingsProviderProps) {
  const defaultSettings: Settings = {
    theme: 'light',
    defaultAtomColor: 'white',
    atomsData: [{ element: 'C', color: 'grey' }]};

  const [settings, setSettings] = useState<Settings>(defaultSettings);

  useEffect(() => {
    const storedSettings = localStorage.getItem('settings');
    if (storedSettings) {
      setSettings(JSON.parse(storedSettings));
    }
  }, []);

  const updateSettings = useCallback((newSettings: Settings) =>{
    setSettings(prevSettings => {
      const updatedSettings = { ...prevSettings, ...newSettings };
      localStorage.setItem('settings', JSON.stringify(updatedSettings));
      return updatedSettings;
    });
  }, []);

  const value = useMemo(() => ({ settings, updateSettings }), [settings, updateSettings]);
  return (
    <SettingsContext.Provider value={ value }>
      {children}
    </SettingsContext.Provider>
  );
}