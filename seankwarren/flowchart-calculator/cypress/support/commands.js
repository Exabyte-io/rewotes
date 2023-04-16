// ***********************************************
// This example commands.js shows you how to
// create various custom commands and overwrite
// existing commands.
//
// For more comprehensive examples of custom
// commands please read more here:
// https://on.cypress.io/custom-commands
// ***********************************************

import 'cypress-drag-drop';

Cypress.Commands.add('dragAndDrop', { prevSubject: 'element' }, (subject, targetSelector) => {
    const dataTransfer = new DataTransfer();
    cy.wrap(subject)
      .trigger('dragstart', { dataTransfer })
      .get(targetSelector)
      .trigger('drop', { dataTransfer })
      .trigger('dragend');
});

Cypress.Commands.add('connectHandles', (sourceSelector, targetSelector) => {
    cy.get(sourceSelector)
      .then(($source) => {
        return cy.get(targetSelector)
          .then(($target) => {
            return {
              sourceId: $source.attr('data-nodeid'),
              targetId: $target.attr('data-nodeid'),
              sourceHandleId: $source.attr('data-handleid'),
              targetHandleId: $target.attr('data-handleid'),
            };
          });
      })
      .then(({ sourceId, targetId, sourceHandleId, targetHandleId }) => {
        console.log('Source ID:', sourceId);
        console.log('Target ID:', targetId);
        console.log('Source Handle ID:', sourceHandleId);
        console.log('Target Handle ID:', targetHandleId);
  

        cy.window().then((win) => {
            win.handleConnect({
                source: sourceId,
                target: targetId,
                sourceHandle: sourceHandleId,
                targetHandle: targetHandleId
            });
        });
      });
  });