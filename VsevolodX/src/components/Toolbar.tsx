import { Button, Dialog, Divider, FileInput, Menu, MenuItem, Navbar, NavbarGroup } from '@blueprintjs/core'
import React, { useState } from 'react'
import styles from './Toolbar.module.scss';
import { useSourceContext } from './SourceContext';
import { useSettings } from './SettingsContext';
import { Popover2 } from '@blueprintjs/popover2';

const Toolbar = () => {
  const { importSource, setSourceName } = useSourceContext();
  const [isDialogOpen, setIsDialogOpen] = useState(false);
  const theme = useSettings().settings.theme;
  const [isEditShown, setEditShown] = useState(false);

  const handleImport = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      importSource(file);
      setSourceName(file.name);
    }
  };

  return (
    <Navbar className={styles.Toolbar}>
      <NavbarGroup>
           <Dialog className={`bp4-${theme}`} isOpen={isDialogOpen} onClose={() => setIsDialogOpen(false)} title="Import file">
              <FileInput text="Choose file..." onInputChange={handleImport} />
            </Dialog>
            
          <Button 
            icon='layout-grid'
            title='Edit'
          />
          <Divider style={{width: '2rem'}}/>
          {/* TODO: Add other Menus and make them look good*/}
          <Popover2 isOpen={isEditShown} position='bottom-right' 
          content={
          <Menu>
            <MenuItem 
              icon='insert' 
              text='Import file' 
              onClick={() => setIsDialogOpen(true)}>  
              
            </MenuItem>
            <MenuItem icon='import' text='Export file' />
          </Menu>
          }>
            <Button icon='document' onClick={() => setEditShown(v => !v)}>File</Button>
          </Popover2>
          
          <Button icon='edit'>Edit</Button>  
          <Button icon='eye-open'>View</Button> 
          <Button icon='help'>Help</Button>  
      </NavbarGroup>
      <NavbarGroup align='right'>
          <Button icon='person' title='account' intent='primary'/>
      </NavbarGroup>
    </Navbar>

  )
}

export default Toolbar