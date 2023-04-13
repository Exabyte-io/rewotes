import React, { useState } from 'react'
import { Button, Dialog, Divider, FileInput, Navbar, NavbarGroup } from '@blueprintjs/core'
import styles from './Toolbar.module.scss';
import { useSourceContext } from '../../context/SourceContext';
import { useSettings } from '../../context/SettingsContext';
import { FileMenu } from '../file_menu/FileMenu';
import { EditMenu } from '../edit_menu/EditMenu';
import { ViewMenu } from '../view_menu/ViewMenu';
import { HelpMenu } from '../help_menu/HelpMenu';
import MyPopover from '../utils/MyPopover';

interface ToolbarProps {
  toggleSettings: () => void,
}

const Toolbar = ({toggleSettings}: ToolbarProps) => {
  const { importSource, setSourceName } = useSourceContext();
  const [isDialogOpen, setIsDialogOpen] = useState(false);
  const theme = useSettings().settings.theme;
  const [isFileShown, setFileShown] = useState(false);
  const [isEditShown, setEditShown] = useState(false);
  const [isViewShown, setViewShown] = useState(false);
  const [isHelpShown, setHelpShown] = useState(false);

  const handleImport = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      importSource(file);
      setSourceName(file.name);
    }
  };

  const handleUndo = () => {
    console.log('You can`t undo what has been done...');
  }

  const handleRedo = () => {
    console.log('You can`t redo it either...');
  }

  return (
    <Navbar className={styles.Toolbar}>
      <NavbarGroup>

        <Dialog className={`bp4-${theme}`} isOpen={isDialogOpen} onClose={() => setIsDialogOpen(false)} title="Import file">
          <FileInput text="Choose file..." onInputChange={handleImport} />
        </Dialog>

        <Button
          icon='layout-grid'
          title='Menu'
          onClick={() => console.log('Clicked on Menu')}
        />
        <Divider style={{ width: '2rem' }} />

        <MyPopover
          isOpen={isFileShown}
          onInteraction={(v) => setFileShown(v)}
          content={<FileMenu openDialog={() => setIsDialogOpen(true)} />}
        >
          <Button icon="document">File</Button>
        </MyPopover>


        <MyPopover
          isOpen={isEditShown}
          onInteraction={(v) => setEditShown(v)}
          content={<EditMenu
              handleUndo={handleUndo}
              handleRedo={handleRedo}
            />
          }
        >
          <Button icon="edit">Edit</Button>
        </MyPopover>

        <MyPopover
          isOpen={isViewShown}
          onInteraction={(v) => setViewShown(v)}
          content={<ViewMenu />}
        >
          <Button icon="eye-open" >View</Button>
        </MyPopover>

        <MyPopover
          isOpen={isHelpShown}
          onInteraction={(v) => setHelpShown(v)}
          // eslint-disable-next-line react/jsx-no-undef
          content={<HelpMenu />}
        >
          <Button icon="help">Help</Button>
        </MyPopover>

      </NavbarGroup>
      <NavbarGroup align='right'>
        <Button icon='cog' title='Settings' intent='none' onClick={() => toggleSettings()} />
        <Button icon='person' title='Account' intent='primary' />
      </NavbarGroup>
    </Navbar>

  )
}

export default Toolbar