export interface ChemicalElement {
    name: string;
    periodicSymbol: string;
    color: string;
    selectedColor: string;
}

export const Symbol2Element:{[periodicSymbol:string]:ChemicalElement} = {
    'H': { name: 'Hydrogen', periodicSymbol: 'H', color: 'cadetblue', selectedColor: 'darkblue' },
    'C': { name: 'Carbon', periodicSymbol: 'C', color: 'charcoal', selectedColor: 'silver' },
    'N': { name: 'Nitrogen', periodicSymbol: 'N', color: 'green', selectedColor: 'lime' },
    'O': { name: 'Oxigen', periodicSymbol: 'O', color: 'blue', selectedColor: 'darkblue' },
};