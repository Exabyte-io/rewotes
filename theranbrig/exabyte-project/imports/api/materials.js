import { Mongo } from 'meteor/mongo';
import { check } from 'meteor/check';

export const Materials = new Mongo.Collection('materials');

console.log(Materials);
Meteor.methods({
  'materials.insert'(name, nodes, connections) {
    console.log(name, nodes, connections);
    check(name, String);
    check(nodes, Array);
    check(connections, Array);
    Materials.insert({
      name,
      nodes,
      connections,
    });
  },

  'materials.remove'(materialId) {
    const material = Material.findOne(materialId);
    if (material) {
      Materials.remove(taskId);
    }
  },

  
});

if (Meteor.isServer) {
  Meteor.publish('materials', function () {
    return Materials.find();
  });
}
