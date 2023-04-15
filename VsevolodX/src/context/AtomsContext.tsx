import { createContext, useCallback, useContext, useEffect, useMemo, useState } from 'react';
import { Vector3 } from 'three';
import parseXYZ from '../actions/parseXYZ';
import { useSourceContext } from './SourceContext';

export type Atom = {
    id: number,
    element: string, //TODO: create enum with all the elements
    position: Vector3; // used for easy connection to ThreeJS
}

interface AtomsContextType {
    atoms: Atom[];
    updateAtoms: (atomId: number, newPosition: Vector3) => void;
    saveAtomsToLocalStorage: () => void;
    loadAtomsFromLocalStorage: () => void;
  }

export const AtomsContext = createContext<AtomsContextType | undefined>(undefined);

export function useAtomsContext(): AtomsContextType {
    const context = useContext(AtomsContext);

    if (!context) {
        throw new Error('useAtoms must be used inside AtomsProvider')
    }
    return context;
}

interface AtomsProps {
    children: React.ReactNode;
}

export function AtomsProvider({ children }: AtomsProps) {
    const [atoms, setAtoms] = useState<Atom[]>([]);
    const {source} = useSourceContext();
    console.log('Atoms Provider: ', atoms);

    useEffect(() => {
        const parsedAtoms = parseXYZ(source);
        if (parsedAtoms.isValid && parsedAtoms.atoms) {
            setAtoms(parsedAtoms.atoms);
        }
    }, [source]);

    const updateAtoms = useCallback((atomId: number, newPosition: Vector3) => {
      setAtoms((prevAtoms) => {
        const updatedAtoms = prevAtoms.map((atom) =>
          atom.id === atomId ? { ...atom, position: newPosition } : atom
        );

        return updatedAtoms;
      });
    }, []);
  
    const saveAtomsToLocalStorage = useCallback(() => {
      localStorage.setItem('Atoms', JSON.stringify(atoms));
    }, [atoms]);
  
    const loadAtomsFromLocalStorage = useCallback(() => {
      const storedAtoms = localStorage.getItem('Atoms');
      if (storedAtoms) {
        setAtoms(JSON.parse(storedAtoms));
      }
    }, []);
  
    const value = useMemo(
      () => ({ atoms, updateAtoms, saveAtomsToLocalStorage, loadAtomsFromLocalStorage }),
      [atoms, updateAtoms, saveAtomsToLocalStorage, loadAtomsFromLocalStorage]
    );
  
    return (
    <AtomsContext.Provider value={value}>
        {children}
    </AtomsContext.Provider>
    );
  };