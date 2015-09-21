import FWCore.ParameterSet.Config as cms

patJetCorrFactorsAK8 = cms.EDProducer("JetCorrFactorsProducer",
                                      src = cms.InputTag("ak8PFJetsCHS"),
                                      emf = cms.bool(False),
                                      extraJPTOffset = cms.string('L1FastJet'),
                                      primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                      levels = cms.vstring('L1FastJet',
                                                           'L2Relative',
                                                           'L3Absolute'),
                                      useNPV = cms.bool(True),
                                      rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                      useRho = cms.bool(True),
                                      payload = cms.string('AK8PFchs'),
                                      flavorType = cms.string('J')
                                      )

patJetsAK8 = cms.EDProducer("PATJetProducer",
                            addJetCharge = cms.bool(False),
                            addGenJetMatch = cms.bool(False),
                            embedGenJetMatch = cms.bool(False),
                            addAssociatedTracks = cms.bool(False),
                            addBTagInfo = cms.bool(False),
                            partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
                            addGenPartonMatch = cms.bool(False),
                            JetPartonMapSource = cms.InputTag(""),
                            resolutions = cms.PSet(),
                            genPartonMatch = cms.InputTag(""),
                            addTagInfos = cms.bool(False),
                            addPartonJetMatch = cms.bool(False),
                            embedGenPartonMatch = cms.bool(False),
                            efficiencies = cms.PSet(),
                            genJetMatch = cms.InputTag(""),
                            useLegacyJetMCFlavour = cms.bool(False),
                            userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
            ),
        userInts = cms.PSet(
            src = cms.VInputTag("")
            ),
        userFloats = cms.PSet(
            src = cms.VInputTag("ak8PFJetsCHSPrunedLinks",
                                "ak8PFJetsCHSSoftDropLinks",
                                "NjettinessAK8:tau1", 
                                "NjettinessAK8:tau2",
                                "NjettinessAK8:tau3")
            ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
            ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring()
        ),
                            jetSource = cms.InputTag("ak8PFJetsCHS"),
                            addEfficiencies = cms.bool(False),
                            discriminatorSources = cms.VInputTag(),
                            trackAssociationSource = cms.InputTag(""),
                            tagInfoSources = cms.VInputTag(cms.InputTag("")),
                            jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsAK8")),
                            embedPFCandidates = cms.bool(False),
                            addJetFlavourInfo = cms.bool(False),
                            addResolutions = cms.bool(False),
                            getJetMCFlavour = cms.bool(False),
                            addDiscriminators = cms.bool(False),
                            jetChargeSource = cms.InputTag(""),
                            JetFlavourInfoSource = cms.InputTag(""),
                            addJetCorrFactors = cms.bool(True),
                            jetIDMap = cms.InputTag(""),
                            addJetID = cms.bool(False)
                            )

selectedPatJetsAK8 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("patJetsAK8"),
    cut = cms.string('pt > 100')
)

redoPatJets = cms.Sequence(patJetCorrFactorsAK8 +
                           patJetsAK8 + 
                           selectedPatJetsAK8)
