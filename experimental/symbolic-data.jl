
using LightXML

function load_system_symbolic_data(filename)
    # parse ex1.xml:
    # xdoc is an instance of XMLDocument, which maintains a tree structure
    xdoc = parse_file(filename)

    # get the root element
    xroot = root(xdoc)  # an instance of XMLElement
    # print its name
    println(name(xroot))  # this should print: bookstore

    # traverse all its child nodes and print element names
    for c in child_nodes(xroot)  # c is an instance of XMLNode
        println(nodetype(c))
        if is_elementnode(c)
            e = XMLElement(c)  # this makes an XMLElement instance
            println(name(e))
        end
    end
end


