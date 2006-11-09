#ifndef __EMPROPERTIES_INTERFACE_H
#define __EMPROPERTIES_INTERFACE_H

class EmPropertiesInterface
{
public:
    EmPropertiesInterface() {}
    virtual ~EmPropertiesInterface() {}

    virtual void read_database(const Dummy& db) = 0;
    virtual void read_database_bowing_parameters(const Dummy& db) = 0;

private:
    EmPropertiesInterface( const EmPropertiesInterface& source );
    void operator = ( const EmPropertiesInterface& source );
};


#endif // __EMPROPERTIES_INTERFACE_H
