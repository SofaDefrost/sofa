/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

#include <string>
namespace sofa::helper
{

class NameDecoder
{
public:
    /// Helper method to get the type name
    template<class T>
    static std::string getTypeName()
    {
        return decodeTypeName(typeid(T));
    }

    /// Helper method to get the class name
    template<class T>
    static std::string getClassName()
    {
        return decodeClassName(typeid(T));
    }

    /// Helper method to get the namespace name
    template<class T>
    static std::string getNamespaceName()
    {
        return decodeNamespaceName(typeid(T));
    }

    /// Helper method to get the template name
    template<class T>
    static std::string getTemplateName()
    {
        return decodeTemplateName(typeid(T));
    }

    template< class T>
    static std::string shortName( const T* = nullptr )
    {
        std::string shortname = getClassName<T>();
        if( !shortname.empty() )
        {
            *shortname.begin() = char(::tolower(*shortname.begin()));
        }
        return shortname;
    }

    /// Helper method to decode the type name
    static std::string decodeFullName(const std::type_info& t);

    /// Helper method to decode the type name to a more readable form if possible
    static std::string decodeTypeName(const std::type_info& t);

    /// Helper method to extract the class name (removing namespaces and templates)
    static std::string decodeClassName(const std::type_info& t);

    /// Helper method to extract the namespace (removing class name and templates)
    static std::string decodeNamespaceName(const std::type_info& t);

    /// Helper method to extract the template name (removing namespaces and class name)
    static std::string decodeTemplateName(const std::type_info& t);
};

}
