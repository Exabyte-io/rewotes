import React from 'react';
import styled from 'styled-components';

const CurrentNodesTable = ({ elements, hoveredElement, removeFromEditor }) => {
  return (
    <TableStyles>
      <h3>Current Points</h3>

      <table className='table-container'>
        <thead>
          <tr>
            <th>ID</th>
            <th>X</th>
            <th>Y</th>
            <th>Z</th>
            <th></th>
          </tr>
        </thead>
        <tbody>
          {elements.map((element, idx) => (
            <tr key={element.id} className={element.id === hoveredElement ? 'hovered' : ''}>
              <td>{element.id}</td>
              <td>{element.x}</td> <td>{element.y}</td>
              <td>{element.z}</td>
              <td>
                <button className='remove-button' onClick={() => removeFromEditor(element.id)}>
                  <img src='/images/delete.svg'/>
                </button>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </TableStyles>
  );
};


export const TableStyles = styled.div`
  td,
  th {
    text-align: center;
    width: 20%;
    padding: 10px 0;
  }
  table {
    margin-top: 20px;
    border-spacing: 0px;
    width: 100%;
  }
  thead {
    display: table-header-group;
    vertical-align: middle;
  }
  th {
    background: #273746;
    font-weight: 300;
  }
  tr:nth-child(even) {
    background: #808b96;
  }
  tr:nth-child(odd) {
    background: #566573;
  }

  tr.hovered {
    background: #2471a3;
    color: white;
  }
  .table-container {
    width: 100%;
    padding: 10px;
    border-radius: 10px;
    background: #273746;
  }
  .remove-button {
    background: transparent;
    box-shadow: none;
    border: none;
    color: white;
  }
  img {
    height: 16px;
  }
`;

export default CurrentNodesTable;
