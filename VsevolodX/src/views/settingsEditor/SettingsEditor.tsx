import React from 'react';
import { useSettings } from '../../context/SettingsContext';
import { Button, Card, Icon, MenuItem } from '@blueprintjs/core';
import styles from './SettingsEditor.module.scss';
import Menu from '../../components/menu/Menu';
import VStack from '../../components/utils/VStack';
import ViewHeading from '../../components/view_heading/ViewHeading';

function SettingsEditor() {
  const { settings, updateSettings } = useSettings();
  const atomsDisplayData = useSettings().settings.atomsDisplayData || [];
  const editingIn3D = useSettings().settings.editingIn3D;

  function toggleTheme() {
    const newTheme = settings.theme === 'light' ? 'dark' : 'light';
    updateSettings({ theme: newTheme });
  }

  async function loadAtomsDisplayData() {
    const res = await fetch('/atomsDisplayData.json');
    const data = await res.json();
    console.log('Atoms display:', data);
    updateSettings({atomsDisplayData: data.atomsDisplayData});
  }

  function toggleEditMode() {
    updateSettings({editingIn3D: !editingIn3D})
  }

  return (
    <Card className={styles.SettingsEditor}>
      <VStack>
      <ViewHeading>
      <h4>Settings</h4>
      </ViewHeading>

      <Button 
        onClick={toggleTheme}
        icon={settings.theme === 'light'? 'moon' : 'lightbulb'}
      >
        {settings.theme === 'light'? 'Dark' : 'Light'} Theme
      </Button>
      {/*TODO: Add ability to change colors in Settings */}
      <Button onClick={() => loadAtomsDisplayData()} >Load atoms preferences</Button>
      <Button onClick={() => toggleEditMode()} >Edit {editingIn3D? 'Source' : '3D'}</Button>
      <Menu>
        {atomsDisplayData.map((element) => {
          return(
            <MenuItem 
            key={element.element}
            text={element.element} 
            icon={<Icon icon='full-circle' color={element.color} />}
            />
            )
          })}
      </Menu>
      </VStack>
    </Card>
  );
}

export default SettingsEditor;
