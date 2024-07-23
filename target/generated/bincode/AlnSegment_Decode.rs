impl :: bincode :: Decode for AlnSegment
{
    fn decode < __D : :: bincode :: de :: Decoder > (decoder : & mut __D)
    ->core :: result :: Result < Self, :: bincode :: error :: DecodeError >
    {
        let variant_index = < u32 as :: bincode :: Decode >:: decode(decoder)
        ?; match variant_index
        {
            0u32 =>Ok(Self ::FullMatch {}), 1u32
            =>Ok(Self ::Match
            {
                0 : :: bincode :: Decode :: decode(decoder) ?, 1 : :: bincode
                :: Decode :: decode(decoder) ?,
            }), 2u32
            =>Ok(Self ::Insertion
            { 0 : :: bincode :: Decode :: decode(decoder) ?, }), variant
            =>Err(:: bincode :: error :: DecodeError :: UnexpectedVariant
            {
                found : variant, type_name : "AlnSegment", allowed : &::
                bincode :: error :: AllowedEnumVariants :: Range
                { min: 0, max: 2 }
            })
        }
    }
} impl < '__de > :: bincode :: BorrowDecode < '__de > for AlnSegment
{
    fn borrow_decode < __D : :: bincode :: de :: BorrowDecoder < '__de > >
    (decoder : & mut __D) ->core :: result :: Result < Self, :: bincode ::
    error :: DecodeError >
    {
        let variant_index = < u32 as :: bincode :: Decode >:: decode(decoder)
        ?; match variant_index
        {
            0u32 =>Ok(Self ::FullMatch {}), 1u32
            =>Ok(Self ::Match
            {
                0 : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?, 1
                : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?,
            }), 2u32
            =>Ok(Self ::Insertion
            { 0 : :: bincode :: BorrowDecode :: borrow_decode(decoder) ?, }),
            variant
            =>Err(:: bincode :: error :: DecodeError :: UnexpectedVariant
            {
                found : variant, type_name : "AlnSegment", allowed : &::
                bincode :: error :: AllowedEnumVariants :: Range
                { min: 0, max: 2 }
            })
        }
    }
}