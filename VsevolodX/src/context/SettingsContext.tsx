import { createContext, useCallback, useContext, useEffect, useMemo, useState } from 'react';

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

  interface SettingsProviderProps {
    children: React.ReactNode;
  }
  

export function SettingsProvider({ children } : SettingsProviderProps) {
  const [settings, setSettings] = useState<Settings>({});

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