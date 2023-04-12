import { useState, useEffect, useMemo, useCallback } from 'react';
import { Settings, SettingsContext } from './SettingsContext';


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
