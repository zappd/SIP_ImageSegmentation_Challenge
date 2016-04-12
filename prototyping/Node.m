classdef Node
   properties
      adj_nodes
   end
   methods
      function obj = Node()
         obj.adj_nodes = Node.empty();
      end
   end
end