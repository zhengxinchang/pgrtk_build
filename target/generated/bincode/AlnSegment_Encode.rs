impl :: bincode :: Encode for AlnSegment
{
    fn encode < __E : :: bincode :: enc :: Encoder >
    (& self, encoder : & mut __E) ->core :: result :: Result < (), :: bincode
    :: error :: EncodeError >
    {
        match self
        {
            Self ::FullMatch
            =>{
                < u32 as :: bincode :: Encode >:: encode(& (0u32), encoder) ?
                ; Ok(())
            }, Self ::Match(field_0, field_1)
            =>{
                < u32 as :: bincode :: Encode >:: encode(& (1u32), encoder) ?
                ; :: bincode :: Encode :: encode(field_0, encoder) ?; ::
                bincode :: Encode :: encode(field_1, encoder) ?; Ok(())
            }, Self ::Insertion(field_0)
            =>{
                < u32 as :: bincode :: Encode >:: encode(& (2u32), encoder) ?
                ; :: bincode :: Encode :: encode(field_0, encoder) ?; Ok(())
            },
        }
    }
}