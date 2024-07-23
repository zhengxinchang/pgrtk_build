impl :: bincode :: Encode for CompactSeq
{
    fn encode < __E : :: bincode :: enc :: Encoder >
    (& self, encoder : & mut __E) ->core :: result :: Result < (), :: bincode
    :: error :: EncodeError >
    {
        :: bincode :: Encode :: encode(&self.source, encoder) ?; :: bincode ::
        Encode :: encode(&self.name, encoder) ?; :: bincode :: Encode ::
        encode(&self.id, encoder) ?; :: bincode :: Encode ::
        encode(&self.seq_frag_range, encoder) ?; :: bincode :: Encode ::
        encode(&self.len, encoder) ?; Ok(())
    }
}