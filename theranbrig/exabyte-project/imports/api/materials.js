import { Mongo } from 'meteor/mongo';
import { check } from 'meteor/check';

export const Materials = new Mongo.Collection('materials');

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
    check(materialId, String);

    const material = Materials.findOne({ _id: materialId });
    console.log(material);
    if (!material) {
      throw new Meteor.error('No Material Found');
    }

    Materials.remove(materialId);
  },

  'materials.update'(materialId, name, nodes, connections) {
    check(materialId, String);
    check(name, String);
    check(nodes, Array);
    check(connections, Array);

    const material = Materials.findOne({ _id: materialId });

    if (!material) {
      throw new Meteor.error('No Material Found');
    }

    Materials.update(materialId, {
      $set: {
        name,
        nodes,
        connections,
      },
    });
  },
});

if (Meteor.isServer) {
  Meteor.publish('materials', function () {
    return Materials.find();
  });
}
