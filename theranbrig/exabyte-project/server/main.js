import { Meteor } from 'meteor/meteor';
import { Materials } from '/imports/api/materials';

function insertMaterial(name, nodes, connections) {
  Materials.insert({ name, nodes, connections });
}

Meteor.startup(() => {
  // If the Links collection is empty, add some data.
  if (Materials.find().count() === 0) {
    insertMaterial('Big Material', [], []);
  }
});
